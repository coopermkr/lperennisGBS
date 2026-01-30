#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis Phylogenetics
#' @date 2026-01-16
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load libraries
library(tidyverse)
library(ape)
library(Biostrings)
library(ggtree)
library(treeio)

# Load in all sample trees
sampInv <- read.tree(file = "4.phylo/NEsamps/noOut.raxml.supportTBE")
sampBoot <- read.tree(file = "4.phylo/NEsamps/noOut.raxml.bootstraps")

# Load in consensus trees of all populations
consInv <- read.tree(file = "4.phylo/withFL/lpPops.raxml.supportTBE")
consBoot <- read.tree(file = "4.phylo/withFL/lpPops.raxml.bootstraps")

# Plot ML trees
ggtree(sampInv) +
  geom_tiplab(size = 2) +
  geom_treescale() +
  ggtitle("Lupinus perennis phylogeny")

# Add highlights around each clade
sampHighlights <- ggtree(sampInv) +
  geom_tiplab(size = 2) +
  geom_treescale() +
  geom_text(aes(label = node), hjust = 0.5)
  geom_hilight(node = 39, fill = "#C200FB") + # Florida
  geom_hilight(node = 42, fill = "#E2A3C7") + # Massachusetts
  geom_hilight(node = 44, fill = "#DC136C") + # Midatlantic
  geom_hilight(node = 52, fill = "#778da9") + # New Hampshire
  geom_hilight(node = 71, fill = "#EC7D10") + # New York
  geom_hilight(node = 71, fill = "#63A46C")   # Vermont

# Assign groups for plotting
groupInfo <- data.frame(label = sampInv$tip.label,
                        group = str_split_i(sampInv$tip.label, pattern = "-", 1)) |>
  mutate(group = str_replace(group, "../1.align/sorted.trim.", ""))

sampGroup <- as.treedata(sampInv) |>
  full_join(groupInfo, by = "label")

pal <- c("black", "#778da9", "red", "#EC7D10", "#63A46C", 
         "blue", "purple", "#E2A3C7", "yellow")

# Plot tree with color assigned to group
sampColors <- ggtree(sampGroup, mapping = aes(color = group),
                     size = 0.75) +
  geom_tiplab(size = 3) +
  geom_treescale() +
  # Supply custom scale matched to PCA
  scale_color_manual(values = pal) +
  labs(title = "Maximum Likelihood L. perennis Phylogeny",
       subtitle = "Estimated with 23,327 SNPs") +
  guides(color = "none",
         size = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 14)) +
  geom_text(aes(label=node))

# Add column showing which nodes to display bootstrap support for
dispNodes <- fortify(sampGroup) |>
  filter(node %in% c(103, 104, 105, 107, 110)) |>
  mutate(boot = round(as.numeric(label), digits = 3))

# Add the root bootstraps (automatically 100%)
rootNodes <- fortify(sampGroup) |>
  filter(node %in% c(103, 104, 105, 107, 110)) |>
  mutate(boot = as.numeric(label))

rootNodes$boot[1] <- 1
rootNodes$boot[2] <- 1

sampColors |>
  # Rearrange nodes
  flip(89, 166)
  rotate(88)
  #flip(104, 103) +
  #rotate(42) |>
  #rotate(43) +
  #geom_text(aes(label=node))
  # Add node labels for important bootstrap numbers
  # geom_text(data = dispNodes,
  #           mapping = aes(label = boot),
  #           color = "steelblue",
  #           #nudge_x = -0.07,
  #           #nudge_y = -4.75,
  #           size = 6) +
  # hexpand(0.05)

jpeg(filename= "4.phylo/salPhylo.jpg", width = 1000, height = 800, quality = 200)
salPhylo
dev.off()
