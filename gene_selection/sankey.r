#devtools::install_github("davidsjoberg/ggsankey")

library(xlsx)
library(ggsankey)
library(ggplot2)
library(dplyr)
library(tidyr)

#setwd("D:/Colabs/Gabriela/Sankey/")

##################################################################
# Import

dir = "~/Documentos/Mestrado/candida_lungs/candida-infection-pipeline"
mod = "turquoise"

load(paste0(dir, "/results/save/gba_lncrna_", mod, ".RData"))


##################################################################
# https://rpubs.com/techanswers88/sankey-with-own-data-in-ggplot
# https://www.youtube.com/watch?v=XRu_Nb8hfIA

# Step 1

df = tabela %>%
  make_long(external_gene_name, go)
df

df$node = factor(df$node,levels = unique(df$node[order(df$node,decreasing = T)]))
df$next_node = factor(df$next_node,levels = unique(df$next_node[order(df$next_node,decreasing = T)]))


# Step 2
dagg = df%>%
  dplyr::group_by(node)%>%
  tally()

# dagg =dagg%>%
#   dplyr::group_by(node)%>%
#   dplyr::mutate(pct = n/TotalCount)


# Step 3
df2 = merge(df, dagg, by.x = 'node', by.y = 'node', all.x = TRUE)
head(df2)

ggplot(df2, aes(x = x,
                next_x = next_x,
                node = node,
                next_node = next_node,
                fill = factor(node),
                label = paste0(node ))) + #,' (', n,')'
  geom_sankey(flow.alpha = .6, space = 1, width = 0.20) +#
  geom_sankey_text(size = 3, color = "black",space = 1,width = 0.1,hjust = 1.5) + #
  scale_fill_viridis_d(option = "viridis",alpha = .5,direction = -1,begin = 0.5,end = 0.9) + #scale_fill_manual
  theme_alluvial(base_size = 18) +
  labs(x = NULL) + theme_bw() + 
  theme(legend.position = "none",axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

