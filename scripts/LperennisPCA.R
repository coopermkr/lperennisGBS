#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis PCA
#' @date 2024-07-31
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load required libraries
library(adegenet)
library(vcfR)
library(tidyverse)

#### Calculate PCA
# Read in the VCF we just saved
missvcf <- read.vcfR("2.stacks/refiltered.allpops.vcf")

# Prepare population info
popList <- read.delim(file = "2.stacks/refmap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2) |>
  filter(id %in% colnames(missvcf@gt))

# Convert VCF to a genInd object and provide population assignments
lpInd <- vcfR2genind(missvcf, pop = popList$pop)
lpInd <- lpInd[indNames(lpInd) != c("../1.align/sorted.trim.92-HK-19", "../1.align/sorted.trim.92-HK-21")] 
lpInd <- lpInd[indNames(lpInd) != c("../1.align/sorted.trim.PM-10", "../1.align/sorted.trim.PM-5", "../1.align/sorted.trim.PM3")]


# Mean missing data
x.lp <- tab(lpInd, freq = TRUE, NA.method = "mean")

# Calculate PCA
pca.lp <- dudi.pca(x.lp, center = TRUE, scale = FALSE, scannf = FALSE, nf = 125)

# Save eigens
eigCoords <- pca.lp$li |>
  cbind(popList) |>
  mutate(region = pop)

write.csv(eigCoords, file = "3.pca/lpEigens.csv")


#### Graph PCA

# First, here are our populations grouped by state from north to south:
# CL, 92-CL
# AB, CN, HK, 92-HK
# SA, AL, MO
# C, P, M, F, PM
# ANF, BW, NK

eigCoords <- read.csv(file = "3.pca/lpEigens.csv")

eig.perc <- 100*pca.lp$eig/sum(pca.lp$eig)

# Plot
palPop <- c("#DC136C","#FFD166", "#F78764", "#63A46C", "#700548",
            "#E2A3C7","#023047","#AD343E","#EC7D10","#33FFFF","#DE9E36",
            "#219ebc","#478978","#C200FB","#778da9",  "#58A4B0","#63D471")

lpPops <- ggplot(data = eigCoords, mapping = aes(x = Axis1, y = Axis2, color = pop)) +
  geom_point(size = 2) +
  scale_color_manual(values = palPop) +
  ggtitle("L. perennis PCA by populations \nwith 20,184 SNPs") +
  labs(title = "L. perennis PCA by population \n with 46,000 SNPs",
       x = "Principal Component 1 (50.6%)",
       y = "Principal Component 2 (3.14%)") +
  guides(size = "none") +
  theme_classic()

#png(filename= "lpAllPops.png", width = 900, height = 500)
#lpPops
#dev.off()

# If we want to simplify and graph by our 6 regions:
eigCoords$region <- gsub("VT1-92", "Vermont", x = eigCoords$region)
eigCoords$region <- gsub("VT1-22", "Vermont", x = eigCoords$region)
eigCoords$region <- gsub("NH3", "New Hampshire", x = eigCoords$region)
eigCoords$region <- gsub("NH1", "New Hampshire", x = eigCoords$region)
#eigCoords$region <- gsub("92-HK", "New Hampshire", x = eigCoords$region)
eigCoords$region <- gsub("NH2-22", "New Hampshire", x = eigCoords$region)
eigCoords$region <- gsub("NY2", "New York", x = eigCoords$region)
eigCoords$region <- gsub("NY1", "New York", x = eigCoords$region)
eigCoords$region <- gsub("MA1", "Massachusetts", x = eigCoords$region)
eigCoords$region <- gsub("FL1", "Florida", x = eigCoords$region)
eigCoords$region <- gsub("FL3", "Florida", x = eigCoords$region)
eigCoords$region <- gsub("FL2", "Florida", x = eigCoords$region)
eigCoords$region <- gsub("MI2", "Midwest/Midatlantic", x = eigCoords$region)
eigCoords$region <- gsub("MD2", "Midwest/Midatlantic", x = eigCoords$region)
eigCoords$region <- gsub("PA1S", "Midwest/Midatlantic", x = eigCoords$region)
eigCoords$region <- gsub("IN2", "Midwest/Midatlantic", x = eigCoords$region)
#eigCoords$region <- gsub("PM", "Midwest/Midatlantic", x = eigCoords$region)

palFl <- c("#C200FB","#E2A3C7", "#DC136C", "#778da9", "#EC7D10", "#63A46C")

lpRegs <- ggplot(data = eigCoords, mapping = aes(x = Axis1, y = Axis2, color = region)) +
  geom_point(size = 3) +
  scale_color_manual(values = palFl) +
  labs(subtitle = "PCA of 40,118 SNPs",
       x = "Principal Component 1 (50.6%)",
       y = "Principal Component 2 (3.14%)",
       tag = "(A)") +
  guides(size = "none", color = guide_legend(title = "Region")) +
  theme_light(base_size = 16) +
  stat_ellipse(type = "norm", linetype = 2) +
  theme(legend.position=c(.725, .81))

png(filename= "3.pca/lpRegions.png", width = 6.5, height = 6.5, units = "in", res = 100)
lpRegs
dev.off()

#### Without Florida ####
lpvcf <- read.vcfR("2.stacks/allPopsFiltered.vcf")

# Convert VCF to a genInd object and provide population assignments
lpInd <- vcfR2genind(lpvcf, pop = popList$pop)

lpInd$pop # We want to keep 1:33, 52:104, 116:136

nofl <- lpInd[c(1:35, 54:113, 125:142),]

# check it worked
nofl$pop

# Mean missing data
x.nofl <- tab(nofl, freq = TRUE, NA.method = "mean")

# Calculate PCA
pca.nofl <- dudi.pca(x.nofl, center = TRUE, scale = FALSE, scannf = FALSE, nf = 100)

eig.perc <- 100*pca.nofl$eig/sum(pca.nofl$eig)

# Save eigens
noflEigs <- data.frame(coords = pca.nofl$li,
                       pops = nofl$pop,
                       Region = nofl$pop)

# If we want to simplify and graph by our 6 regions:
noflEigs$Region <- gsub("92-CL", "Vermont", x = noflEigs$Region)
noflEigs$Region <- gsub("CL", "Vermont", x = noflEigs$Region)
noflEigs$Region <- gsub("AB", "New Hampshire", x = noflEigs$Region)
noflEigs$Region <- gsub("CN", "New Hampshire", x = noflEigs$Region)
noflEigs$Region <- gsub("92-HK", "New Hampshire", x = noflEigs$Region)
noflEigs$Region <- gsub("HK", "New Hampshire", x = noflEigs$Region)
noflEigs$Region <- gsub("SA", "New York", x = noflEigs$Region)
noflEigs$Region <- gsub("AL", "New York", x = noflEigs$Region)
noflEigs$Region <- gsub("MO", "Massachusetts", x = noflEigs$Region)
noflEigs$Region <- gsub("ANF", "Florida", x = noflEigs$Region)
noflEigs$Region <- gsub("NK", "Florida", x = noflEigs$Region)
noflEigs$Region <- gsub("BW", "Florida", x = noflEigs$Region)
noflEigs$Region <- gsub("^C$", "Midwest/Midatlantic", x = noflEigs$Region)
noflEigs$Region <- gsub("^PM$", "Midwest/Midatlantic", x = noflEigs$Region)
noflEigs$Region <- gsub("^F$", "Midwest/Midatlantic", x = noflEigs$Region)
noflEigs$Region <- gsub("^M$", "Midwest/Midatlantic", x = noflEigs$Region)
noflEigs$Region <- gsub("^P$", "Midwest/Midatlantic", x = noflEigs$Region)

palReg <- c("#E2A3C7", "#DC136C", "#778da9", "#EC7D10", "#63A46C")

# Plot PCA
noflPCA <- ggplot(data = noflEigs, mapping = aes(x = coords.Axis1, y = coords.Axis2, color = Region)) +
  geom_point(size = 3) +
  scale_color_manual(values = palReg) +
  labs(subtitle = "PCA of 20,184 SNPs",
       x = "Principal Component 1 (6.84%)",
       y = "Principal Component 2 (6.50%)") +
  guides(size = "none", color = guide_legend(title = "Region")) +
  theme_light(base_size = 16) +
  stat_ellipse(type = "norm", linetype = 2) +
  theme(legend.position=c(.8,.2))

png(filename= "3.pca/noflPCA.png", width = 6.5, height = 6.5, units = "in", res = 100)
noflPCA
dev.off()

#write.csv(noflEigs, file = "3.pca/noflEigs.csv")

#### Just Northeast SNPs ####
nevcf <- read.vcfR("2.stacks/neFiltered.vcf")

# Prepare population info
neList <- read.delim(file = "2.stacks/nemap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2)

# Convert VCF to a genInd object and provide population assignments
neInd <- vcfR2genind(nevcf, pop = neList$pop)

neInd <- lpInd[row.names(lpInd@tab) %in% neList$id]

# Mean missing data
x.ne <- tab(neInd, freq = TRUE, NA.method = "mean")

# Calculate PCA
pca.ne <- dudi.pca(x.ne, center = TRUE, scale = FALSE, scannf = FALSE, nf = 83)

eig.perc <- 100*pca.ne$eig/sum(pca.ne$eig)

# Save eigens
neEigs <- data.frame(coords = pca.ne$li,
                     pops = neList$pop,
                     Region = neList$pop,
                     year = neList$pop) |>
  mutate(Region = str_remove_all(pops, '[123456789-]'),
         year = as.factor(case_when(year == "VT1-92" ~ 92, .default = 22)))

# If we want to simplify and graph by our 6 regions:
neEigs$Region <- gsub("92-CL", "Vermont", x = neEigs$Region)
neEigs$Region <- gsub("CL", "Vermont", x = neEigs$Region)
neEigs$Region <- gsub("AB", "New Hampshire", x = neEigs$Region)
neEigs$Region <- gsub("CN", "New Hampshire", x = neEigs$Region)
neEigs$Region <- gsub("92-HK", "New Hampshire", x = neEigs$Region)
neEigs$Region <- gsub("HK", "New Hampshire", x = neEigs$Region)
neEigs$Region <- gsub("SA", "New York", x = neEigs$Region)
neEigs$Region <- gsub("AL", "New York", x = neEigs$Region)
neEigs$Region <- gsub("MO", "Massachusetts", x = neEigs$Region)

# And make a column to add shapes for collection years:
neEigs$year <- gsub("92-CL", "1992", x = neEigs$year)
neEigs$year <- gsub("CL", "2022", x = neEigs$year)
neEigs$year <- gsub("AB", "2022", x = neEigs$year)
neEigs$year <- gsub("CN", "2022", x = neEigs$year)
neEigs$year <- gsub("HK", "2022", x = neEigs$year)
neEigs$year <- gsub("SA", "2022", x = neEigs$year)
neEigs$year <- gsub("AL", "2022", x = neEigs$year)
neEigs$year <- gsub("MO", "2022", x = neEigs$year)

neEigs <- read_csv("3.pca/neEigs.csv")

palReg <- c("#E2A3C7", "#778da9", "#EC7D10", "#63A46C")

# Plot PCA
nePCA <- ggplot(data = neEigs, mapping = aes(x = coords.Axis1, y = coords.Axis2, 
                                             color = Region, shape = year)) +
  geom_point(size = 3) +
  scale_shape_manual(values = c(16, 17)) +
  scale_color_manual(values = palReg) +
  labs(subtitle = "PCA of 10,945 SNPs",
       x = "Principal Component 1 (10.6%)",
       y = "Principal Component 2 (6.51%)",
       tag = "(B)") +
  guides(shape = "none", color = guide_legend(title = "Region"), tab = "(B)") +
  theme_light(base_size = 12) +
  stat_ellipse(type = "norm", linetype = 2) +
  theme(legend.position=c(.85,.64))

png(filename= "3.pca/nePCA.png", width = 315, height = 380)
nePCA
dev.off()

#write.csv(neEigs, file = "3.pca/neEigs.csv")

#### Combined Figure ####
library(gridExtra)
library(grid)

jpeg(filename= "3.pca/pcaFig.jpg", width = 11, height = 6.75, units = "in", res = 100)
grid.arrange(lpRegs, nePCA, ncol = 2)
dev.off()


#### Previous Code ####

# Load in popfile from STACKS run
popList <- read.csv(file = "3.pca/noOut.csv")

# Generate a genind file for adegenet
lupine <- vcfR2genind(lpvcf, pop = popList, return.alleles = TRUE)

# Check for missing data
sum(is.na(lupine$tab))
# There are X NAs

# Replace NAs with mean allele frequency
lupine_data <- scaleGen(lupine, NA.method = "mean")
# Check it worked
dim(lupine_data)

save(lupine_data, file = "lupineDataMean.Rdata")

# Generate PCA
lupinePCA1 <- dudi.pca(lupine,
                       cent = FALSE, scale = FALSE, SCANNF = FALSE)

# Calculate Axis variance
vals <- read.delim(file = "3.pca/noOutEvals.ssv", sep = " ")
variance <- 100*vals/sum(vals)
head(variance)
# Axis1: 4.37%
# Axis2: 3.14%
# Axis3: 2.75%
# Axis4: 2.65%

# Combine the eigenvectors with population info
eigen <- read.delim(file = "3.pca/noOutEvecs.ssv", header= TRUE, sep = " ") |>
  cbind(popList) |>
  select(!X) |>
  rename(pop = V1)

#### Vizualize
# First, here are our populations grouped by state from north to south:
# CL, 92-CL
# AB, CN, HK, 92-HK
# SA, AL, MO
# C, P, M, F, PM
# ANF, BW, NK

# And here's the order they output in (alphabetical):
unique(popList$V1)

# Now here are the matching colors:

pal <- c("#DC136C","#FFD166", "#F78764", "#63A46C", "#700548",
  "#E2A3C7","#023047","#AD343E","#EC7D10","#33FFFF","#DE9E36",
  "#219ebc","#478978","#C200FB","#778da9",  "#58A4B0","#63D471")
         
# Plot, taking out that one crazy outlier point
#Everything: Super clear distinction with lat. on Axis1 and long. on Axis2
noOutAll <- ggplot(data = eigen, mapping = aes(x = Axis1, y = Axis2, color = pop)) +
  geom_point(size = 2) +
  scale_color_manual(values = pal) +
  ggtitle("L. perennis PCA by populations \nwith 1.4 million SNPs") +
  guides(size = "none") +
  theme_classic()

png(filename= "noOutAll.png", width = 900, height = 500)
noOutAll
dev.off()


# Axes 1 and 2 by population sets
#Florida
ggplot(data= eigen |> filter(V1 %in% c("AL", "NK", "BW", "ANF")), 
       mapping = aes(x= Axis1, y = Axis2, color = V1)) +
  geom_point()
  

#Old and New
ggplot(data= eigen |> filter(V1 %in% c("CL", "92-CL", "AL", "HK", "92-HK", "AB")), 
       mapping = aes(x= Axis1, y = Axis2, color = V1)) +
  geom_point()

#Midwest
ggplot(data= eigen |> filter(V1 %in% c("AL", "C", "F", "M", "P", "PM")), 
       mapping = aes(x= Axis1, y = Axis2, color = V1)) +
  geom_point()

#Northeast
palNE <- c("#63A46C", "#AD343E", "#EC7D10", "#478978", "#63D471")

ggplot(data= eigen |> filter(V1 %in% c("AL", "CL", "CN", "MO", "SA")), 
       mapping = aes(x= Axis1, y = Axis2, color = V1)) +
  geom_point() +
  geom_point(size = 2) +
  scale_color_manual(values = palNE) +
  ggtitle("L. perennis PCA by populations \nwith 1.4 million SNPs") +
  guides(size = "none")

# Axes 3 and 4 by population sets
#Everything: little pattern
ggplot(data = eigen, mapping = aes(x = Axis3, y = Axis4, color = V1)) +
  geom_point()

#Florida
ggplot(data= eigen |> filter(V1 %in% c("AL", "NK", "BW", "ANF")), 
       mapping = aes(x= Axis3, y = Axis4, color = V1)) +
  geom_point()

#Old and new
ggplot(data= eigen |> filter(V1 %in% c("CL", "92-CL", "AL", "HK", "92-HK", "AB")), 
       mapping = aes(x= Axis3, y = Axis4, color = V1)) +
  geom_point()

#Midwest
ggplot(data= eigen |> filter(V1 %in% c("AL", "C", "F", "M", "P", "PM")), 
       mapping = aes(x= Axis3, y = Axis4, color = V1)) +
  geom_point()

#Northeast (repeat CL)
ggplot(data= eigen |> filter(V1 %in% c("AL", "CN", "SA", "MO", "CL")), 
       mapping = aes(x= Axis3, y = Axis4, color = V1)) +
  geom_point()




