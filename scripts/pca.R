#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis PCA Northeast Sites
#' @date 2024-07-31
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load required libraries
library(adegenet)
library(vcfR)
library(tidyverse)
library(cowplot)

#### PCA for Northeast Sites
nevcf <- read.vcfR("2.stacks/neFiltered.vcf")

neList <- read.delim(file = "2.stacks/nemap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2)

# Convert to genind
neind <- vcfR2genind(nevcf, pop = neList$pop, return.alleles = TRUE)

# Pull out tab
netab <- tab(neind, freq=TRUE, NA.method = "mean")

# Calculate PCA
nePCA <- dudi.pca(netab,
                       center = TRUE, scale = TRUE,
                       scannf = FALSE, nf = 200)

vec <- nePCA$li
val <- nePCA$eig
write.table(x= vec, file = "3.pca/neEvecs.ssv")
write.table(x = val, file = "3.pca/neEvals.ssv")

# Calculate explained variance
var <- data.frame(Axis = 1:77,
                  Variance = 100*val/sum(val))

head(var) # Taking out that crazy outlier might help this!
write.table(var, file = "axisVariance.ssv")

# Make a bar chart of the variance explained
varplot <- ggplot(data = var,
       mapping = aes(x = Axis, y = Variance)) +
  geom_col() +
  theme_classic() +
  labs(y = "% Variance Explained")

# Plot the PCA eigenvectors
vec <- vec |>
  mutate(region = substr(rownames(vec), 24, 25),
         region = str_replace_all(region, c("92"="Seed Bank", "AB"="NH", "AL"="NY", "CL"="VT", "CN"="NH", "HK"="NH", "MO"="MA", "SA"="NY")))
vec$region[1:6] <- "VT"
vec$region[7:8] <- "NH"


palMap <- c("#E2A3C7", "#778da9", "#EC7D10", "#63A46C")

nepca <- ggplot(data = vec,
       mapping = aes(x = -Axis1, y = Axis2, color = region)) +
  geom_point(size = 2) +
  theme_classic(base_size = 16) +
  theme(plot.title = element_text(hjust = 0.5)) +
  labs(x = "Axis 1 (5.34%)",
    y = "Axis 2 (4.41%)") +
  scale_color_manual(values = palMap, name = "State")

jpeg(filename = "3.pca/nepca.jpg", height = 6, width = 7, units = "in", res = 300)
nepca
dev.off()
