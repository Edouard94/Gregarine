# Load libraries
library(ggtree)
library(ggplot2)
library(ape)
library(tidytree)
library(treeio)
library(dplyr)
library(extrafont)
loadfonts()
setwd("C:\\Users\\eb826\\OneDrive - University of Exeter\\Documents")

# Read in your tree file (replace "treefile.nwk" with the path to your tree file)
tree <- read.tree("Greg18S-Tree-Clade1.tree")

# Define the outgroup
outgroup <- "FJ459737-Amoebogregarina-nigra"
# Reroot the tree with "Chromera_velia" as the outgroup
tree <- root(tree, outgroup)

# Define a minimum for the branch length
#branch_length_min <- 0.001

# This shortens your tree to fit tip labels. Adjust the factor for a better fit.
xlim_adj <- max(ggtree(tree)$data$x) * 2

# Extend the length of your branches by multiplying the edge lengths by a factor (e.g., 1.5)
tree$edge.length <- tree$edge.length * 1

# Plot the tree with new labels
p <- ggtree(tree, ladderize = TRUE) + 
  geom_tiplab(hjust = 0, size = 5, linesize=.5, offset=0.001, fontface="italic", family="Times New Roman") + 
  geom_treescale(y = -0.95, fontsize = 3.9) +
  geom_text2(aes(label = round(as.numeric(label), 2), 
                 subset = !is.na(as.numeric(label)) & as.numeric(label) > 0 & as.numeric(label) <= 1), 
             vjust = -0.5, hjust = -0.2, size = 4, check_overlap = TRUE) + 
  theme(legend.text = element_text(size=8)) + 
  xlim(0, xlim_adj) +
  scale_fill_identity(guide = FALSE)

# Display the tree
p

ggsave("Greg18S-Tree-Clade1.svg", p, dpi = 600, width = 7, height = 5, units = "in")
