#dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
dir = "~/workspace/candida_lungs/candida-infection-pipeline"

setwd(dir)

options(stringsAsFactors = FALSE)

df = read.csv(file = paste0(dir, "/results/MMxGS_expression.csv"))

tabela_ml = read.table(file.path(paste0(dir, "/results/27K-ranking.csv")), header = FALSE, sep = ",")
genes_ml = subset(tabela_ml, tabela_ml$V2 == 1)

genes = subset(df, df$ensembl_gene_id %in% genes_ml$V1 &
                abs(df$log2FoldChange) > 1 & df$padj < 0.05 &
                (df$GS.condition < -0.6 | df$GS.condition > 0.6) & df$MM > 0.6)

modules = c("turquoise", "blue")
edges_file = NULL
new_edges_file = NULL
nodes_file = NULL
ext_info = NULL

filter_edges_file = function(){
  for(module in modules){
    edges_file = read.table(file = paste0(dir, "/results/CytoscapeInput-edges-", module, ".txt"),
                                      header = TRUE, sep = "\t")
    
    lncrna = subset(genes, genes$gene_biotype == "lncrna" & genes$moduleColor == module)
    not_lncrna = lncrna = subset(genes, genes$gene_biotype != "lncrna" & genes$moduleColor == module)
    selected_genes = subset(genes, genes$moduleColor == module)

    new_edges_file = subset(edges_file, 
                       (edges_file$fromNode %in% selected_genes$ensembl_gene_id &
			edges_file$toNode %in% selected_genes$ensembl_gene_id))    

    #new_edges_file = subset(edges_file, 
    #                   (edges_file$fromNode %in% lncrna$ensembl_gene_id & edges_file$toNode %in% lncrna$ensembl_gene_id) |
    #                    (edges_file$fromNode %in% lncrna$ensembl_gene_id & edges_file$toNode %in% not_lncrna$ensembl_gene_id) |
    #                      (edges_file$fromNode %in% not_lncrna$ensembl_gene_id & edges_file$toNode %in% lncrna$ensembl_gene_id))
    
    write.table(new_edges_file, file = paste0(dir, "/results/CytoscapeInput-edges-", module, "-filtered.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

add_external_info = function(){
  for(module in modules){
    nodes_file = read.table(file = paste0(dir, "/results/CytoscapeInput-nodes-", module, ".txt"),
                                          header = TRUE, sep = "\t")
    temp = subset(df, df$ensembl_gene_id %in% nodes_file$nodeName)
    info = c("ensembl_gene_id", "gene_biotype", "GS.condition", "MM")
    temp = temp[info]
    ext_info = merge(nodes_file, temp, by = 1)
    
    write.table(ext_info, file = paste0(dir, "/results/CytoscapeInput-nodes-", module, "-ext.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}

filter_edges_file()
#add_external_info()
