library(ggplot2)

repository = "candida-infection-pipeline"

# For RStudio only
current_dir = dirname(rstudioapi::getActiveDocumentContext()$path)
dir = strsplit(current_dir, repository)[[1]][1]

res_total = read.csv(file = paste0(dir, repository, "/results/MMxGS_expression.csv"))
res = subset(res_total, res_total$moduleColor %in% c("turquoise", "blue"))
rownames(res) = res$ensembl_gene_id

biotypes = NULL
biotypes[["PROT"]] = c("protein_coding", "Protein coding")
biotypes[["LNCRNA"]] = c("lncRNA", "Long non-coding RNA")

biotype = biotypes[["LNCRNA"]]

res = subset(res, res$gene_biotype == biotype[1])

sig_t = rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1 & res$moduleColor == "turquoise"))
sig_t_color = rep("turquoise", length(sig_t))
names(sig_t_color) = sig_t

sig_b = rownames(subset(res, padj < 0.05 & abs(log2FoldChange) > 1 & res$moduleColor == "blue"))
sig_b_color = rep("blue", length(sig_b))
names(sig_b_color) = sig_b

nsig = rownames(subset(res, !rownames(res) %in% c(sig_t, sig_b)))
nsig_color = rep("grey60", length(nsig))
names(nsig_color) = nsig

colScale = c(sig_t_color, nsig_color, sig_b_color)
colScale = colScale[order(match(names(colScale), rownames(res)))]

res_plot = subset(res, res$ensembl_gene_id %in% c(sig_t, nsig, sig_b))
y_plot = -log10(res_plot$padj)
padj_max = max(y_plot[is.finite(y_plot)])

ggplot(data = res_plot, aes(x = log2FoldChange, y = -log10(padj), col = colScale, label = NA)) + 
    geom_point(shape = 19, alpha = 0.6, size = 1.5) +
    theme_minimal() +
    xlim(-8, 8) +  ylim(0, padj_max) + 
    ggtitle("Volcano plot", subtitle = paste0(biotype[2], " genes")) +
    xlab(expression(Log[2]*" fold change")) + ylab(expression(-Log[10]* " adjusted p-value")) + 
    scale_x_continuous(breaks = scales::pretty_breaks(n = 7)) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
    theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
          plot.title = element_text(size = 16, hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 10),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 10)) + 
    scale_color_manual(name = "Differential expression", values = c("blue", "grey70", "turquoise"),
                       labels = c("Turquoise module DEGs", "Non-significant", "Blue module DEGs")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "grey45") +
    geom_vline(xintercept = -1, linetype = "dashed", colour = "grey45") +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey45") 

ggsave(
    filename = paste0("volcano_plot_mod_", biotype[1], ".pdf"),
    plot = last_plot(),
    path = paste0(dir, repository, "/results/figures"),
    scale = 1,
    width = 6.5,
    height = 6,
    units = "in",
    dpi = 300,
    bg = "white"
)
