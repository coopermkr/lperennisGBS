#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis PCA
#' @date 2024-07-31
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

setwd("/project/pi_brook_moyers_umb_edu/lupine/")
# Load required libraries
library(adegenet)
library(vcfR)
library(dplyr)
library(ggplot2)

#### Process VCF for visualization
lpvcf <- read.vcfR("3.pca/noOut.snps.vcf")

# Make a popfile (note this is edited beforehand to not include BW1-2)
#popList <- as.data.frame(read.delim("2.stacks/pops", header = FALSE))
#write.csv(popList, file = "3.pca/noOut.csv")
pops <- read.csv("3.pca/noOut.csv", header = TRUE) |> select(V1)

# Generate a genind file for adegenet
lupine <- vcfR2genind(lpvcf, pop = pops, return.alleles = TRUE)

# Check for missing data
sum(is.na(lupine$tab))
# There are 315984422 NAs (?) I feel like that's wrong because that's basically the whole file

# Replace NAs with mean allele frequency
#lupineMean <- scaleGen(lupine, NA.method = "mean")
# Check it worked
#class(lupineMean)
#dim(lupineMean) # 140 samples, 2799556 rows

# Let PCA handle missing data
x.lupine <- tab(lupine, freq=TRUE, NA.method = "mean")
write.table(x.lupine, file = "3.pca/meanNoOutLupine.genind")

#### Calculate PCA
lupinePCA1 <- dudi.pca(x.lupine,
                       center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = 200)
# Find clusters
#clust <- find.clusters(lupine, max.n.clust=20)

# Calculate DAPC
#dapc1 <- dapc(lupine, clust$grp)

# Save eigenvectors and values
vec <- lupinePCA1$li
val <- lupinePCA1$eig
write.table(x= vec, file = "3.pca/noOutEvecs.ssv")
write.table(x = val, file = "3.pca/noOutEvals.ssv")

# Calculate explained variance
var <- 100*lupinePCA1$eig/sum(lupinePCA1$eig)

head(var) # Taking out that crazy outlier might help this!
write.table(var, file = "axisVariance.ssv")

# Also I'm not exactly mad given how much data there is here
lupinePCA1$eig[1]

#### Visualize
#s.label(lupinePCA1$li)
#lPCA <- s.class(lupinePCA1$li, fac = pop(lupine), col=funky(17))
# The integrated functions for viz don't work, so we'll make our own with the eigenvectors

save(lupinePCA1, file = "3.pca/noOutPCA1.R")
