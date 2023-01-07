library(WGCNA)

options(stringsAsFactors = FALSE)
options(digits = 15)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
#dir = "~/workspace/candida_lungs/candida-infection-pipeline"

lnames = load(file = paste0(dir, "/results/save/mm-dataInput.RData"))
lnames

topology_analysis = function(){
	# Choose a set of soft-thresholding powers
	powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
	# network topology analysis function
	sft = pickSoftThreshold(datExpr, powerVector = powers, networkType = "signed", verbose = 5, blockSize = ncol(datExpr))
	
	sizeGrWindow(5, 5)
	par(mfrow = c(1, 2))
	cex1 = 0.9
	
	# Scale-free topology fit index as a function of the soft-thresholding power
	png(file = paste0(dir, "/results/network_topology.png"))
	plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit,signed R^2", type = "n",
	     main = paste("Scale independence"))

	text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
	     labels = powers, cex = cex1, col = "red")

	abline(h = 0.90, col = "red")
	dev.off()

	png(file = paste0(dir, "/results/network_mean_connectivity.png"))
	plot(sft$fitIndices[,1], sft$fitIndices[,5],
	     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
		main = paste("Mean connectivity"))
	text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = cex1, col = "red")
	dev.off()

	save(sft, file = paste0(dir, "/results/save/mm-topologyAnalysis.RData"))	

	return(sft)
}

#sft = topology_analysis()

lnames = load(file = paste0(dir, "/results/save/mm-topologyAnalysis.RData"))
lnames

power_value = 20

# One-step network construction and module detection
net = blockwiseModules(datExpr = datExpr, power = power_value, maxBlockSize = ncol(datExpr),
                       networkType = "signed",
                       TOMType = "signed", minModuleSize = 50,
                       reassignThreshold = 0, mergeCutHeight = 0.3,
                       numericLabels = FALSE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = paste0("mmTOM-p", power_value),
                       verbose = 3)
table(net$colors)

sizeGrWindow(12, 9)
pdf(file = paste0(dir, "/results/figures/network_modules_p", power_value, ".pdf"),
    width = 10,
    height = 4.5)
plotDendroAndColors(net$dendrograms[[1]], net$colors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, 
                    autoColorHeight = FALSE, colorHeight = 0.25)
dev.off()

moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]

save(sft, net, MEs, moduleColors, geneTree,
     file = paste0(dir, "/results/save/mm-networkConstruction-p", power_value, ".RData"))
