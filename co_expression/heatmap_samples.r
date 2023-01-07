library(WGCNA)
library(stringr)

options(stringsAsFactors = FALSE)
options(digits = 15)

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
#dir = "~/workspace/candida_lungs/candida-infection-pipeline"

lnames = load(file = paste0(dir, "/results/save/mm-dataInput.RData"))
lnames
lnames = load(file = paste0(dir, "/results/save/mm-networkConstruction-p20.RData"))
lnames

rede_mod = read.csv(file = paste0(dir, "/results/MMxGS_expression.csv"))

matrizexpressao = read.csv(file = paste0(dir, "/results/expression_matrix.csv"))

matrizexpressao = t(subset(matrizexpressao, matrizexpressao$X %in% rede_mod$ensembl_gene_id))
colnames(matrizexpressao) = matrizexpressao[1, ] 

matrizexpressao = as.data.frame(matrizexpressao[-c(1), ])
matrizexpressao = as.data.frame(sapply(matrizexpressao, as.numeric))

sampleTable = read.table(file.path(dir, "metadata/SraRunTable.txt"), header = TRUE, sep = ",")

#__________________________________________
datME = MEs
signif(cor(datME, use = "p"), 2)

par(mar = rep(2, 4))
plotMEpairs(datME)

table(sampleTable$Condition)

Control = rownames(datME)[c(1:11)]
Candida = rownames(datME)[c(12:19)]

#__________________________________________

datTraits$Condition[datTraits$Condition == 0] = "Control"
datTraits$Condition[datTraits$Condition == 1] = "Disease"

traits = model.matrix(~0+Condition, datTraits)
traits = as.data.frame(traits)

moduleTraitCor = cor(MEs, traits, use = "p")

#__________________________________________

# Candida

colorList = unique(moduleColors)
labelList = c(paste0("CTRL", 1:11), paste0("CAND", 1:8))

sizeGrWindow(7, 8)
pdf(file = paste0(dir, "/results/modules_vs_samples.pdf"),
    width = 6,
    height = 6)  

for(color in colorList){
    ME = datME[, paste("ME", color, sep = "")]
    
    par(mfrow = c(2, 1), mar = c(0.3, 5.5, 4, 2))
    plotMat(t(scale(matrizexpressao[ , net$colors == color])),
            nrgcols = 10, rlabels = F, rcols = color, #clabels = labelList,
            main = paste0("Gene expression vs. \n", str_to_title(color),
                          " module eigengene expression"),
            cex.main = 1) 
    
    colors = data.frame(Samples = rownames(datME))
    
    mod = paste0("ME", color)
    
    if(moduleTraitCor[mod, "ConditionControl"] > moduleTraitCor[mod, "ConditionDisease"]){
        colors$colors = ifelse(colors$Samples %in% Candida, color, "gray87")
    }else{
        colors$colors = ifelse(colors$Samples %in% Candida, "gray87", color) #"orchid"
    }
    
    par(mar = c(5, 4.7, 0, 1.2)) #c(bottom, left, top, right)
    barplot(ME, col = colors$colors, main = "", cex.main = 1,
            ylab = "Eigengene expression", xlab = "Samples",
            names.arg = labelList, las = 2, cex.names = 0.7, space = 0.15)
}

dev.off()
