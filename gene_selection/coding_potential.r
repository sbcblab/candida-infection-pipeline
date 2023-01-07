library(readr)
library(tibble)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"

modulos = c("turquoise", "blue")
exp = NULL
cpc = NULL
res_cpc = NULL
res_exp_cpc = NULL

for(mod in modulos){
    exp[[mod]] = read.csv(file = paste0(dir, "/results/selected_genes_", mod, ".csv"), sep = "\t")
    #cpc[[mod]] = read.csv(file = paste0(dir, "/results/nova/result_cpc2_", mod, ".txt"), sep = "\t")
    cpc[[mod]] = read.csv(file = paste0(dir, "/results/nova/cpc2_cdna_", mod, ".txt"), sep = "\t")
    
    gene = NULL
    transcript = NULL
    gene_name = NULL
    
    for(x in 1:nrow(cpc[[mod]])){
        temp = unlist(strsplit(cpc[[mod]][x, "X.ID"], "_"))
        gene = append(gene, temp[1])
        transcript = append(transcript, temp[2])
        temp2 = subset(exp[[mod]], exp[[mod]]$ensembl_gene_id == temp[1])
        gene_name = append(gene_name, temp2$external_gene_name)
    }
    
    cpc[[mod]] = add_column(cpc[[mod]], gene, .after = 1)
    cpc[[mod]] = add_column(cpc[[mod]], transcript, .after = 2)
    cpc[[mod]] = add_column(cpc[[mod]], gene_name, .after = 3)
    cpc[[mod]]$X.ID = NULL
    
    coding = subset(cpc[[mod]], cpc[[mod]]$label == "coding")
    res_cpc[[mod]] = subset(cpc[[mod]], cpc[[mod]]$gene %in% coding$gene)
    res_exp_cpc[[mod]] = subset(exp[[mod]],
                                exp[[mod]]$ensembl_gene_id %in% res_cpc[[mod]]$gene)
}

res_samba = NULL
res_exp_samba = NULL
samba = NULL

for(mod in modulos){
    samba[[mod]] = read.csv(file = paste0(dir, "/results/nova/rnasamba_cdna_", mod, ".tsv"), sep = "\t")
    
    gene = NULL
    transcript = NULL
    gene_name = NULL
    
    for(x in 1:nrow(samba[[mod]])){
        temp = unlist(strsplit(samba[[mod]][x, "sequence_name"], "_"))
        gene = append(gene, temp[1])
        temp2 = unlist(strsplit(temp[2], " "))
        transcript = append(transcript, temp2[1])
        gene_name = append(gene_name, temp2[2])
    }
    
    samba[[mod]] = add_column(samba[[mod]], gene, .after = 1)
    samba[[mod]] = add_column(samba[[mod]], transcript, .after = 2)
    samba[[mod]] = add_column(samba[[mod]], gene_name, .after = 3)
    samba[[mod]]$sequence_name = NULL
    
    coding = subset(samba[[mod]], samba[[mod]]$classification == "coding")
    res_samba[[mod]] = subset(samba[[mod]], samba[[mod]]$gene %in% coding$gene)
    res_exp_samba[[mod]] = subset(exp[[mod]],
                                exp[[mod]]$ensembl_gene_id %in% res_samba[[mod]]$gene)
}

rm(gene, transcript, gene_name, temp, temp2, x, coding, mod)

mod = "blue"
write.table(cpc[[mod]], file = paste0(dir, "/results/nova/cpc2_cdna_", mod, ".tab"),
        sep = "\t", quote = FALSE, row.names = FALSE)



