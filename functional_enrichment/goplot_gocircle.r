#install.packages('GOplot')
library(GOplot)
library(readr)
library(ggplot2)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
expression = read.csv(file = paste0(dir, "/results/MMxGS_expression.csv"))
ontology = read.csv(file = paste0(dir, "/results/topgo_selected_genes_final.tab"), sep = "\t")

#lnames = load(file = paste0(dir, "/results/save/goplot.RData"))

modules = c("turquoise", "blue")
genes = NULL
genes_mod = NULL

for(mod in modules){
    temp = read.csv(file = paste0(dir, "/results/selected_genes_0.4_", substr(mod, 1, 1), ".tsv"), sep = "\t")
    genes_mod[[mod]] = temp$ensembl_gene_id
}

genes = append(genes_mod[[modules[1]]], genes_mod[[modules[2]]])

selected_genes = subset(expression, expression$ensembl_gene_id %in% genes)
expr = data.frame(ID = selected_genes$ensembl_gene_id, gene = selected_genes$external_gene_name,
                       logFC = selected_genes$log2FoldChange)

circ = circle_dat(ontology, expr)
# Reduce redundant terms with a gene overlap >= 0.75...
# reduced_circ = reduce_overlap(circ, overlap = 0.75)

processes = "inflammatory response
chemokine-mediated signaling pathway
positive regulation of tumor necrosis factor production
neutrophil chemotaxis
positive regulation of angiogenesis
positive regulation of ERK1 and ERK2 cascade
extracellular matrix organization
antimicrobial humoral immune response mediated by antimicrobial peptide
negative regulation of apoptotic process
killing of cells of another organism
response to molecule of fungal origin
response to wounding"
processes = unlist(strsplit(processes, "\n"))

pdf(file = paste0(dir, "/results/figures/gocircle.pdf"),
    width = 14,
    height = 8)

GOCircle(circ, nsub = processes, zsc.col = c("red", "white", "blue"),
         rad2 = 3)
dev.off()

#save(ontology, expression, genes, genes_mod, file = paste0(dir, "/results/save/goplot.RData"))
