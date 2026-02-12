#'''''''''''''''''''''''''''''''''
#' Allele Frequency Processing
#' @date 2026-01-21
#' @author Cooper Kimball-Rhines
#'''''''''''''''''''''''''''''''''

library(tidyverse)

# Read in and join all of the allele frequencies with new population column
afConcat <- function(pop){
  afPop <- read_tsv(paste("5.translocation/allelefreqs/", pop, sep = "")) |>
    mutate(population = paste(pop))
}

popList <- list.files("5.translocation/allelefreqs/")

afAll <- map(.x = popList, .f = afConcat) |> list_rbind()

# Split minor and major alleles and filter out NAs
afSplit <- afAll |>
  mutate(chromPos = paste(CHROM, POS, sep = "-")) |> # Make a unique ID
  mutate(major = str_split_i(`{ALLELE:FREQ}`, "\t", 1),
         minor = str_split_i(`{ALLELE:FREQ}`, "\t", 2),
         major = as.numeric(str_split_i(major, ":", 2)),
         minor = as.numeric(str_split_i(minor, ":", 2)),
         population = str_split_i(population, ".pop.", 1)) |> #Split the major and minor AF
  filter(!is.na(minor),
         !is.na(major)) # Filter out NAs

# Check that all loci appear in all populations after filtering
missSites <- afSplit |>
  select(chromPos, `{ALLELE:FREQ}`) |>
  group_by(chromPos) |> count() |> # Count how many populations show up at each locus
  filter(n != 9) # Filter loci that aren't present for all sites
# For some reason these 20 appear twice, so they should be filtered out.

# Find variant loci between populations
varSites <- afSplit |>
  select(chromPos, `{ALLELE:FREQ}`) |>
  group_by(chromPos) |> unique() |> count() |> # Count how many rows are unique based on that ID
  filter(n > 1) # 1 = invariant site, Max is 9 = highly variant site

# Filter the original list to only include variant sites that appear in all populations
afFilt <- afSplit |>
  filter(chromPos %in% varSites$chromPos,
         !chromPos %in% missSites$chromPos)

afFilt |> filter(chromPos == "Scaffold_10__1_contigs__length_23337219-10003245")
# Thankfully the order that the alleles appear is alphabetized, not sorted by major/minor
# So the allele frequencies will be consistent as long as the file is sorted properly

# total number of loci: 25760
afFilt |> group_by(population) |> count()

# Write out frequency lists for each population pairing of interest
# Note: These files must be edited to add the total number of variants and alleles per locus at the top
# recipient population first, donor second (or subsequent)
### Extant VT population
afFilt |> ungroup() |>
  filter(population == "VT22") |>
  arrange(population, chromPos) |>
  select(major, minor) |>
  write_delim(file = "5.translocation/allelefreqs/VT.af", delim = " ", col_names = FALSE)

# VT and seed bank
afFilt |> ungroup() |>
  filter(population %in% c("VT22", "VT92")) |>
  arrange(population, chromPos) |>
  select(major, minor) |>
  write_delim(file = "5.translocation/allelefreqs/VTtransloc.af", delim = " ", col_names = FALSE)

### All NH populations supplemented by AL
Auburn <- afFilt |> ungroup() |>
    filter(population == "AB") |>
    arrange(population, chromPos)

Hooksett <- afFilt |> ungroup() |>
  filter(population == "HK22") |>
  arrange(population, chromPos)

Concord <- afFilt |> ungroup() |>
  filter(population == "CN") |>
  arrange(population, chromPos)

Albany <- afFilt |> ungroup() |>
  filter(population == "AL") |>
  arrange(population, chromPos)

Saratoga <- afFilt |> ungroup() |>
  filter(population == "SA") |>
  arrange(population, chromPos)

Montague <- afFilt |> ungroup() |>
  filter(population == "MA") |>
  arrange(population, chromPos)

# For Concord/Albany real life translocation:
CNAL <- rbind(Concord, Albany) |>
  select(major, minor)

CNAL |>
  write_delim(file = "5.translocation/allelefreqs/CNAL.af", delim = " ", col_names = FALSE)

# Bind together: Three supplemented pops, three solo pops, donor pop
# Albany supplement
NHAL <- rbind(Auburn, Hooksett, Concord, Auburn, Hooksett, Concord, Albany)
# Should be 25760 * 7 = 180320 rows

NHAL |> select(major, minor) |>
  write_delim(file = "5.translocation/allelefreqs/NHAL.af", delim = " ", col_names = FALSE)

## All following simulations only need to run the supplemented pops since the solos count for everyone
# Saratoga supplement
NHSA <- rbind(Auburn, Hooksett, Concord, Saratoga)

NHSA |> select(major, minor) |>
  write_delim(file = "5.translocation/allelefreqs/NHSA.af", delim = " ", col_names = FALSE)

# Montague supplement
NHMO <- rbind(Auburn, Hooksett, Concord, Montague)

NHMO |> select(major, minor) |>
  write_delim(file = "5.translocation/allelefreqs/NHMO.af", delim = " ", col_names = FALSE)
