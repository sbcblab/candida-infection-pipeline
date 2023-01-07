library(WGCNA)
library("RColorBrewer")
library(pals)

options(stringsAsFactors = FALSE)
options(digits = 15)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
#dir = "~/workspace/candida_lungs/candida-infection-pipeline"

lnames = load(file = paste0(dir, "/results/save/mm-dataInput.RData"))
lnames
lnames = load(file = paste0(dir, "/results/save/mm-networkConstruction-p20.RData"))
lnames

datTraits$Condition[datTraits$Condition == 0] = "Control"
datTraits$Condition[datTraits$Condition == 1] = "Disease"

traits = model.matrix(~0+Condition, datTraits)
traits = as.data.frame(traits)

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

moduleColors = net$colors
MEs = net$MEs
geneTree = net$dendrograms[[1]]

# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, traits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

cores = c("green4", "white", "orchid4")
paleta = colorRampPalette(cores)

sizeGrWindow(10, 6)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(10, 10, 3, 3))

pdf(file = paste0(dir, "/results/figures/module_condition.pdf"),
    width = 11,
    height = 3.5)

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = t(moduleTraitCor),
               xLabels = names(MEs), #names(traits),
               xSymbols = names(MEs),
               yLabels = c("Control\ncondition", "Infection\ncondition"),
               ySymbols = c("Control\ncondition", "Infection\ncondition"),
               colorLabels = FALSE,
               colors = paleta(50), #blueWhiteRed(50),
               #colors = rev(brewer.piyg(50)),
               textMatrix = t(textMatrix),
               setStdMargins = TRUE,
               cex.text = 1,
               cex.lab = 1,
               textAdj = 0.5,
               zlim = c(-1, 1),
               xLabelsAngle = 45,
               main = paste("Module-trait Relationships"))
dev.off()

# Gene Significance and Module Membership

condition = as.data.frame(traits$ConditionControl)
names(condition) = "condition"

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep = "")
names(MMPvalue) = paste("p.MM", modNames, sep = "")

geneTraitSignificance = as.data.frame(cor(datExpr, condition, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(condition), sep = "")
names(GSPvalue) = paste("p.GS.", names(condition), sep = "")


# 3.c Intramodular analysis: identifying genes with high GS and MM

library(stringr)

selectModules = c("turquoise", "blue")

for (module in selectModules){
  column = match(module, modNames)
  moduleGenes = moduleColors == module
  
  sizeGrWindow(10, 10)
  par(mfrow = c(1, 1), mar = c(4, 5.5, 3, 2), family = "Helvetica")
  pdf(file = paste0(dir, "/results/figures/mm_gs_", module, ".pdf"),
      width = 5.5,
      height = 4.5)
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = "Module membership",
                     ylab = "Gene significance for health condition",
                     main = paste("Module Membership vs. Gene Significance \n", str_to_title(module), "module \n"),
                     cex.main = 1, cex.lab = 1, cex.axis = 1, col = module)
  abline(v = 0.6, col = "black", lty = 3, lwd = 1)
  abline(h = 0.6, col = "black", lty = 3, lwd = 1)
  dev.off()
}

names(datExpr)
names(datExpr)[moduleColors == "turquoise"]

annot = read.csv(file.path(dir, "/results/gene_annotation.csv"), stringsAsFactors = FALSE)
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$ensembl_gene_id)

sum(is.na(probes2annot))

# Create the starting data frame
geneInfo0 = data.frame(ensembl_gene_id = probes,
                       external_gene_name = annot$external_gene_name[probes2annot],
                       description = annot$description[probes2annot],
                       entrezgene_id = annot$entrezgene_id[probes2annot],
                       gene_biotype = annot$gene_biotype[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

modOrder = order(-abs(cor(MEs, condition, use = "p")))

# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]])
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep = ""),
                       paste("p.MM.", modNames[modOrder[mod]], sep = ""))
}

geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.condition))
geneInfo = geneInfo0[geneOrder, ]

#write.csv(as.data.frame(geneInfo), row.names = FALSE, file = paste0(dir, "/results/gene_modules_p20.csv"))


selectModules = c("turquoise", "blue")
FiltergenesListC1 = list()

deseq = read.csv(file = paste0(dir, "/results/deseq_results.csv"))

for (k in names(geneTraitSignificance)){
  for (module in selectModules){
    #module = selectModules[4]
    column = match(module, modNames)
    modGenes = moduleColors == module
    GSxMM = data.frame(KO = rownames(geneModuleMembership)[modGenes],
                       ModuleMembership = geneModuleMembership[modGenes, column],
                       GeneSignificance = geneTraitSignificance[modGenes, 1])
    GSxMM = GSxMM[order(GSxMM$ModuleMembership, GSxMM$GeneSignificance, decreasing = T),]
    GSxMM = subset(GSxMM, GSxMM$ModuleMembership > 0.8 & abs(GSxMM$GeneSignificance) > 0.8)
    x = paste(k, module, sep = " vs. ")
    #FiltergenesListC1[[x]] = GSxMM
    FiltergenesListC1[[x]] = merge(GSxMM, deseq, by = 0)
  }
}

str(FiltergenesListC1)
names(FiltergenesListC1)

mmgs = merge(geneInfo0, deseq, by = 1, all.x = TRUE, all.y = FALSE)

mod_membership = c()
p_mod_membership = c()

for(i in 1:nrow(geneInfo0)) { # MM e p.MM apenas do módulo principal
  row = mmgs[i,]
  mod_membership = c(mod_membership, row[[paste0("MM.", row$moduleColor)]])
  p_mod_membership = c(p_mod_membership, row[[paste0("p.MM.", row$moduleColor)]])
}

mmgs0 = subset(mmgs, select = c(1:8))
mmgs0["MM"] = mod_membership
mmgs0["p.MM"] = p_mod_membership

mmgs_deseq = merge(mmgs0, deseq, by = 1, all.x = TRUE, all.y = FALSE)

write.csv(as.data.frame(mmgs), row.names = FALSE, file = paste0(dir, "/results/MMxGS_full.csv"))
write.csv(as.data.frame(mmgs_deseq), row.names = FALSE, file = paste0(dir, "/results/MMxGS_expression.csv"))
