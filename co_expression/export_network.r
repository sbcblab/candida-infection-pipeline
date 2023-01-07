library(WGCNA)

workingDir = "~/workspace/candida_lungs/candida-infection-pipeline"
setwd(workingDir)

options(stringsAsFactors = FALSE)

collectGarbage()

lnames = load(file = paste0(workingDir, "/results/save/mm-dataInput.RData"), verbose = TRUE)
lnames = load(file = paste0(workingDir, "/results/save/mm-networkConstruction-p20.RData"))
lnames = load(file = paste0(workingDir, "/results/save/mmTOM-p20-block.1.RData"), verbose = TRUE)

# Recalculate topological overlap if needed
#TOM = TOMsimilarityFromExpr(datExpr, power = 20, networkType = "signed", TOMType = "signed")
#save(TOM, file = paste0(workingDir, "/data/save/mm-tom.RData"))

print(paste0("TOM - dimensão: ", dim(TOM)))
TOM = as.matrix(TOM)
print(paste0("TOM - dimensão: ", dim(TOM)))

collectGarbage()

# Read in the annotation file
annot = read.csv(file = paste0(workingDir, "/data/gene_annotation.csv"))

# Select modules
modules = c("turquoise", "blue")

for(module in modules){
	# Select module probes
	print(paste0("Començando módulo ", module))
	probes = names(datExpr)
	inModule = is.finite(match(moduleColors, module))
	modProbes = probes[inModule]
	modGenes = annot$external_gene_name[match(modProbes, annot$ensembl_gene_id)]

	# Select the corresponding Topological Overlap
	modTOM = TOM[inModule, inModule]
	dimnames(modTOM) = list(modProbes, modProbes)

	# Export the network into edge and node list files Cytoscape can read
	cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", module, ".txt", sep = ""),
                               nodeFile = paste("CytoscapeInput-nodes-", module, ".txt", sep = ""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule])
	
	print(paste0("Terminando módulo ", module))
}

