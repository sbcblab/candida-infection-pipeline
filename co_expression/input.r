library(WGCNA)

options(stringsAsFactors = FALSE)
options(digits = 15)

#dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
dir = "~/workspace/candida_lungs/candida-infection-pipeline"

mmData = read.csv(file.path(dir, "/results/expression_matrix.csv"), stringsAsFactors = FALSE)
mmData = subset(mmData, X != "Condition")
dim(mmData)
names(mmData)

datExpr0 = as.data.frame(t(mmData[, -c(1:1)]))
datExpr0 = as.data.frame(sapply(datExpr0, as.numeric)) #TESTE
names(datExpr0) = mmData$X
rownames(datExpr0) = names(mmData)[-c(1:1)]

gsg = goodSamplesGenes(datExpr0, verbose = 3) # Columns = genes; rows = samples
gsg$allOK

if (!gsg$allOK){
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Genes removed:", paste(dim((datExpr0)[!gsg$goodGenes])[2])))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Samples removed:", paste(dim((datExpr0)[!gsg$goodSamples])[2])))
  
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

#cCOND1 = gsg$goodGenes
#aCOND1 = adjacency(t(cCOND1),type = "distance")
#kCOND1 = as.numeric(apply(aCOND1, 2 ,sum))-1
#ZCOND1.k = scale(kCOND1)
#thresholdZ.k = -2.5
#COND1outlierColor = ifelse(ZCOND1.k<thresholdZ.k,"red","black")
#sampletreeCOND1 = flashClust(as.dist(1-aCOND1), method = "average")
#COND1datColors = data.frame(outlierC=COND1outlierColor)

#png(file = paste0(workingDir, "/outlier.png"))
#plotDendroAndColors(sampletreeCOND1,groupLabels = names(COND1datColors), colors=COND1datColors, 
#                    main = "Dendrogram and outlier detection")
#dev.off()

#save(sampletreeCOND1, COND1datColors,
#     file = paste0(workingDir, "/co_expression/outlier-sample.RData"))

sampleTree = hclust(dist(datExpr0), method = "average")

sizeGrWindow(12, 9)
par(cex = 0.6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub = "", xlab = "", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

abline(h = 15, col = "red")
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)

keepSamples = (clust == 0)
datExpr = datExpr0 [keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

traitData = read.table(file.path(dir, "metadata/SraRunTable.txt"), header = TRUE, sep = ",")
datTraits = data.frame(Condition = traitData$Condition)
rownames(datTraits) = paste0("sample", 1:19)

datTraits$Condition[datTraits$Condition == "Control"] = 0
datTraits$Condition[datTraits$Condition == "Disease"] = 1

save(datExpr, datTraits, file = paste0(dir, "/results/save/mm-dataInput.RData"))
