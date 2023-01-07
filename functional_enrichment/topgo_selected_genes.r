library(topGO)
library(org.Mm.eg.db)
library(xlsx)
library(WGCNA)

workingDir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
setwd(workingDir)

options(stringsAsFactors = FALSE)

lnames = load(file = paste0(workingDir, "/results/save/mm-networkConstruction-p20.RData"))
lnames

table(net$colors)

allGeneIDs = names(moduleColors) # nome dos genes da rede
#namesModulecolor = unname(moduleColors) # nome dos m?dulos da rede
namesModulecolor = c("turquoise", "blue")

genes = list()
deg_up = list()
deg_down = list()
GenebyModules = list()


genes_t = read.csv(file = paste0(workingDir, "/results/selected_genes_turquoise.csv"), sep = "\t")
list_t = read.csv(file = paste0(workingDir, "/results/selected_genes_0.4_t.tsv"), sep = "\t")
genes_t = subset(genes_t, genes_t$ensembl_gene_id %in% list_t$ensembl_gene_id)

genes_b = read.csv(file = paste0(workingDir, "/results/selected_genes_blue.csv"), sep = "\t")
list_b = read.csv(file = paste0(workingDir, "/results/selected_genes_0.4_b.tsv"), sep = "\t")
genes_b = subset(genes_b, genes_b$ensembl_gene_id %in% list_b$ensembl_gene_id)

genes = rbind(genes_t, genes_b)
rownames(genes) = genes[, 1]

deg_up = subset(genes, genes$log2FoldChange > 1 & genes$pvalue < 0.05)
deg_down = subset(genes, genes$log2FoldChange < -1 & genes$pvalue < 0.05)
# Separadando os genes presentes em cada modulo da condicao

GenesbyModule = (net$colors %in% namesModulecolor)
#Name = paste(h)
tmp = list(allGeneIDs[GenesbyModule])
GenebyModules = tmp


str(GenebyModules) # lista com os genes todos separados de acordo com o modulo
rm(Name, tmp, GenesbyModule)

#Selecionado os genes de interesse e rodando o topGO

geneUniverse = allGeneIDs
OntologyResults = list()


# apenas genes entre os selecionados
genesOfInterest = unlist(genes$ensembl_gene_id) #unlist(GenebyModules[[i]])
geneList = factor(as.integer(geneUniverse %in% genesOfInterest))
names(geneList) = geneUniverse # aqui e onde identificamos o que deve ser analisado em comparacao ao universo

#Analisando os dados com o topgo
myGOdata = new("topGOdata", description = genes$moduleColor,  ontology = "BP", allGenes = geneList,  
               annot = annFUN.org, nodeSize = 5, mapping = "org.Mm.eg.db", ID = "Ensembl")

resultTopgo = list()
resultTopgo = runTest(myGOdata, algorithm = "weight01", statistic = "fisher")
mysummaryTopgo = summary(attributes(resultTopgo)$score <= 0.05)
numsignifT = as.integer(mysummaryTopgo[[3]])
allRes = GenTable(myGOdata, topgoFisher = resultTopgo, numChar = 300, ranksOf = "topgoFisher", topNodes = numsignifT)
allRes$FDR = p.adjust(allRes$topgoFisher, method = "fdr")

#criando uma lista dos GO.ID + genes de interesse que estao presentes no GO.ID
myterms = allRes$GO.ID
mygenes = genesInTerm(myGOdata, myterms)
GenesinTerms = NULL

for (j in 1:length(myterms)){
    myterm = myterms[j]
    mygenesforterm = mygenes[myterm][[1]] [mygenes[myterm][[1]] %in% genesOfInterest == TRUE] # colocar apenas os genes que estao no modulo analisado 
    myNumgenes = length(mygenesforterm)
    mygenesforterm3 = paste(mygenesforterm, collapse = ',')
    
    hubs = rownames(subset(genes, genes$MM > 0.8 & abs(genes$GS.condition) > 0.8)) ### df_kME obtido com a fun??o signedKME do WGCNA
    sighubs = intersect(mygenesforterm, hubs)
    sighubs2 = ifelse(length(sighubs) == 0, NA, paste(sighubs, collapse = ','))
    myNumhubs = length(sighubs)
    
    # # Exp vs Ctrl
    DEGSonto1 = intersect(deg_up$ensembl_gene_id, intersect(mygenes[[j]], sigGenes(myGOdata))) # resultados de expressao diferencial ja filtrados por fdr e lof2fc => aqui s? os UP
    DEGSonto1a = ifelse(length(DEGSonto1) == 0, NA, paste(DEGSonto1, collapse = ','))
    myNumbDEGsonto1 = length(DEGSonto1)
    
    DEGSonto2 = intersect(deg_down$ensembl_gene_id,intersect(mygenes[[j]], sigGenes(myGOdata))) # resultados de expressao diferencial ja filtrados por fdr e lof2fc => aqui s? os DOWN
    DEGSonto2a = ifelse(length(DEGSonto2) == 0, NA, paste(DEGSonto2, collapse = ','))
    myNumbDEGsonto2 = length(DEGSonto2)
    
    GenesinTerms = rbind(GenesinTerms, data.frame(GO.ID = myterm, MyGenes = mygenesforterm3, NumbGenes = myNumgenes, 
                                                  SigHubs = sighubs2, NumbHubs = myNumhubs,
                                                  DEGS_UP = DEGSonto1a, NumbDEGS_UP = myNumbDEGsonto1,
                                                  DEGS_DOWN = DEGSonto2a, NumbDEGS_DOWN = myNumbDEGsonto2))
}

FinalTable = merge(allRes, GenesinTerms, by = "GO.ID")

# Salvando o resultados
#Name = paste0("Modules")
tmp = FinalTable
OntologyResults = tmp


str(OntologyResults)
names(OntologyResults)

# Convertendo os resultados em arquivo de excel (salva resultados para todos os modulos, mas usando o objeto OntologyResults da pra selecionar apenas os modulos de interesse )
Ontofile = createWorkbook(type = 'xlsx')

for (k in 1:length(names(OntologyResults))) {
    sheet1 = createSheet(Ontofile, sheetName = paste(names(OntologyResults)[k]))
    addDataFrame(OntologyResults[[k]], sheet1, row.names = F )
}

saveWorkbook(Ontofile, file = paste0(workingDir, "/results/ontology_selected_genes_hubs.xlsx"))

ontology = NULL

df = as.data.frame(OntologyResults)
pref = "Modules"
new_names = gsub(pref, "", names(df))
names(df) = new_names
ontology = data.frame(Category = rep("BP", dim(df)[1]), 
                             ID = df$GO.ID, Term = df$Term, 
                             Genes = df$MyGenes, adj_pval = df$FDR)

write.table(ontology, row.names = FALSE, sep = "\t" , quote = FALSE,
            file = paste0(workingDir, "/results/topgo_selected_genes_final.tab"))


save(ontology, file = paste0(workingDir, "/results/save/ontology_selected_genes_final.RData"))
