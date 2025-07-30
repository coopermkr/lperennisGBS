#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#'
#' L. perennis VCF Processing
#' @date 2024-07-31
#' @author Cooper Kimball-Rhines
#' 
#'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

# Load required libraries
library(SNPfiltR)
library(vcfR)
library(dplyr)

# Load in VCF
lpvcf <- read.vcfR("2.stacks/populations.snps.vcf")

# Load in popMap
popList <- read.delim(file = "2.stacks/refmap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2)

## Missingness Filtering
# Hard filter to a min depth of 5 and quality of 30
afvcf <- hard_filter(vcfR=lpvcf, depth = 3, gq = 30) |>

# Remove loci with > 2 alleles
  filter_biallelic() |>

# Filter out non-heterozygous loci
 filter_allele_balance()

# Filter out very high coverage loci
#lpvcf <- max_depth(lpvcf) # Create histogram of depths to choose a cutoff
# I don't have a high depth problem, so I'm not filtering for this

# Set missingness cutoff
missing_by_sample(vcfR=afvcf, popmap = popList)
sampvcf <- missing_by_sample(afvcf, cutoff = 0.95) |>

# Remove invariant sites
  min_mac(min.mac = 1)

# Set missingness cutoff by SNP
missing_by_snp(sampvcf)
missvcf <- missing_by_snp(sampvcf, cutoff = 0.8) # There's more stats to do later to confirm this is good
# I'm setting this at 0.8 to keep 20k SNPs for 114 samples

# Write filtered vcf
write.vcf(missvcf, "2.stacks/allPopsFiltered.vcf")

# Write linkage filtered vcf
thinvcf <- distance_thin(missvcf, min.distance = 500)

write.vcf(thinvcf, "2.stacks/allPopsLinkageFiltered.vcf")


#### Northeast Sites Filtering ####
lpvcf <- read.vcfR("2.stacks/populations.snps.vcf")

# Load in popMap
neList <- read.delim(file = "2.stacks/nemap.txt", header = FALSE) |>
  rename(id = V1,
         pop = V2)

# Subset VCF to only Northeast samples
nevcf <- lpvcf[,c("FORMAT",neList$id)]

## Missingness Filtering
# Hard filter to a min depth of 5 and quality of 30
afne <- hard_filter(vcfR=nevcf, depth = 3, gq = 30) |>
  
  # Remove invariant sites from subsetting
  min_mac(min.mac = 1) |>
  
  # Remove loci with > 2 alleles
  filter_biallelic() |>
  
  # Filter out non-heterozygous loci
  filter_allele_balance()

# Set missingness cutoff
missing_by_sample(vcfR=afne, popmap = neList)
sampne <- missing_by_sample(afne, cutoff = 0.95) |>
  
  # Remove invariant sites from subsetting
  min_mac(min.mac = 1)

# Set missingness cutoff by SNP
missing_by_snp(sampne)
missne <- missing_by_snp(sampne, cutoff = 0.75) # There's more stats to do later to confirm this is good
# I'm setting this at 0.8 to keep 10k SNPs for 78 samples

# Write filtered vcf
write.vcf(missne, "2.stacks/neFiltered.vcf")

# Write linkage filtered vcf
thinne <- distance_thin(missne, min.distance = 500)

write.vcf(thinne, "2.stacks/neLinkageFiltered.vcf")

