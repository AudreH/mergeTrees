rm(list = ls())

library(dendextend)
library(ggraph)
library(gridExtra)
library(tidyverse)
library(tidygraph)
library(gtable)
library(grid)

load("hc_list_methods.RData")

gr <- as_tbl_graph(hc_list_methods$MT)
layout_MC <- create_layout(gr, layout = "dendrogram", circular=FALSE, height = height)
layout_MC$Celltype = factor(gsub(".*_","", layout_MC$label), levels = c("CD14", "CD4", "CD8", "WB"))
layout_MC$node_color = as.numeric(layout_MC$Celltype)+1
layout_MC$node_size = ifelse(layout_MC$leaf, 4, 0)
layout_MC$label = ""

gMC = ggraph::ggraph(layout_MC) +
  geom_node_point(colour = layout_MC$node_color,
                  size= layout_MC$node_size)+
    geom_node_text(aes(label = label)) +
  geom_edge_elbow() +
  theme_classic() +
  ggtitle("Merged Clustering") +
  xlab("") +
  ylab("Height")+
  theme(axis.text.x=element_blank(), axis.line.x=element_blank(),
        axis.ticks.x = element_blank()) 
gMC


gr <- as_tbl_graph(hc_list_methods$AD)
layout_AD <- create_layout(gr, layout = "dendrogram", circular=FALSE, height = height)
layout_AD$Celltype = factor(gsub(".*_","", layout_AD$label), levels = c("CD14", "CD4", "CD8", "WB"))
layout_AD$node_color = as.numeric(layout_AD$Celltype) + 1
layout_AD$node_size = ifelse(layout_AD$leaf, 4, 0)
layout_AD$label = ""

gAD = ggraph::ggraph(layout_AD) +
  geom_node_point(colour = layout_AD$node_color,
                  size= layout_AD$node_size)+
  geom_node_text(aes(label = label)) +
  geom_edge_elbow() +
  theme_classic() +
  ggtitle("Averaged Distance") +
  xlab("") +
  ylab("Height")+
  theme(axis.text.x=element_blank(), axis.line.x=element_blank(),
        axis.ticks.x = element_blank()) 
gAD

hc_list_methods$DC$labels = hc_list_methods$AD$labels
gr <- as_tbl_graph(hc_list_methods$DC)
layout_DC <- create_layout(gr, layout = "dendrogram", circular=FALSE, height = height)
layout_DC$Celltype = factor(gsub(".*_","", layout_DC$label), levels = c("CD14", "CD4", "CD8", "WB"))
layout_DC$node_color = as.numeric(layout_DC$Celltype) +1
layout_DC$node_size = ifelse(layout_DC$leaf, 4, 0)
layout_DC$label = ""

gDC = ggraph::ggraph(layout_DC) +
  geom_node_point(colour = layout_DC$node_color,
                  size= layout_DC$node_size)+
  geom_node_text(aes(label = label)) +
  geom_edge_elbow() +
  theme_classic() +
  ggtitle("Direct Clustering") +
  xlab("") +
  ylab("Height")+
  theme(axis.text.x=element_blank(), axis.line.x=element_blank(),
        axis.ticks.x = element_blank()) 
gDC

grid.arrange(gMC, gAD, gDC, ncol = 3)

DC_grob <- ggplotGrob(gDC)
AD_grob <- ggplotGrob(gAD)
MC_grob <- ggplotGrob(gMC)

all_grob <- gtable_cbind(DC_grob, AD_grob)
all_grob <- gtable_cbind(all_grob, MC_grob)

grid.newpage()
grid.draw(all_grob)

ggsave(all_grob, filename = 'Figure4.pdf', width = 30, height = 10, units = 'cm')

